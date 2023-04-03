import requests, auth

# replace with your Pinata API key and secret
PINATA_API_KEY = "40d47b677e0c7cec6456"
PINATA_SECRET_API_KEY = "2af04f84c6e380691844f7eebc765be533270540a09e20abaa4d20128bc512f3"

# make API request to get list of pinned files
response = requests.get('https://api.pinata.cloud/data/pinList?pageLimit=100', 
                        headers={'pinata_api_key': PINATA_API_KEY, 'pinata_secret_api_key': PINATA_SECRET_API_KEY})

# check if API request was successful
if response.status_code == 200:
    # extract CIDs from response data
    pinned_files = response.json()['rows']
    cids_to_unpin = [file['ipfs_pin_hash'] for file in pinned_files]
    
    # iterate over CIDs and unpin each file
    for cid in cids_to_unpin:
        response = requests.delete(f'https://api.pinata.cloud/pinning/unpin/{cid}',
                                headers={'pinata_api_key': PINATA_API_KEY, 'pinata_secret_api_key': PINATA_SECRET_API_KEY})
        if response.status_code == 200:
            print(f'{cid} was unpinned successfully')
        else:
            print(f'Error unpinning {cid}: {response.text}')
else:
    print('Error fetching pinned files:', response.text)
